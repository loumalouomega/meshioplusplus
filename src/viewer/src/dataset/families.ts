/**
 * The training-spec model families, as the launch form sees them.
 *
 * A pure twin of `physicsnemo/_train.py`'s `_FAMILIES` table: which spec
 * block each family reads, which form panel shows its options, and which
 * `train_start` keyword arguments the panel's inputs feed. The panel posts
 * ONLY the chosen family's keys — the server refuses another family's by
 * name, so posting a UI default for every family would turn a default into
 * an error the user never asked for.
 *
 * The family `name`s here are pinned equal, in order, to the Python
 * `_MODELS` tuple by a test on the Python side (`tests/python/
 * test_train_spec.py`), the `manifest.ts` parity precedent: a family added
 * on one side and not the other is a red build, not a silently missing
 * option.
 */

export type FamilyBlock = 'Graph' | 'Grid' | 'Operator';

/** How a form control's text becomes the tool argument. */
export type FieldKind =
    /** a plain number; blank means "leave the default" */
    | 'number'
    /** three positive integers, `16,16,16` — a grid resolution */
    | 'triple'
    /** two positive integers, `8,8` — an AFNO patch */
    | 'pair'
    /** a world axis 0/1/2 from a select; blank means "no squeeze" */
    | 'axis'
    /** comma-separated non-empty names (the `Metadata` keys of a DeepONet) */
    | 'names'
    /** a select's value, verbatim */
    | 'choice';

export interface FamilyField {
    /** the `train_start` keyword argument */
    readonly key: string;
    /** the form control's element id */
    readonly id: string;
    readonly kind: FieldKind;
    /** blank is an error rather than "omit" */
    readonly required?: boolean;
}

export interface Family {
    readonly name: string;
    readonly label: string;
    readonly block: FamilyBlock;
    /** the `<div>` holding this family's options */
    readonly panelId: string;
    readonly fields: readonly FamilyField[];
    /** the `train_defaults.frameworks` key that must be true to offer it,
     * beyond `physicsnemo` itself */
    readonly requires?: string;
}

export const FAMILIES: readonly Family[] = [
    {
        name: 'meshgraphnet',
        label: 'MeshGraphNet (mesh graph)',
        block: 'Graph',
        panelId: 't-graph-opts',
        fields: [
            { key: 'processor_size', id: 't-processor', kind: 'number' },
            { key: 'hidden_dim', id: 't-hidden', kind: 'number' },
        ],
    },
    {
        name: 'srresnet',
        label: 'SRResNet (superresolution grid)',
        block: 'Grid',
        panelId: 't-grid-opts',
        fields: [
            { key: 'resolution', id: 't-resolution', kind: 'triple', required: true },
            { key: 'scaling_factor', id: 't-scaling', kind: 'number' },
            { key: 'conv_layer_size', id: 't-conv', kind: 'number' },
            { key: 'resid_blocks', id: 't-blocks', kind: 'number' },
        ],
    },
    {
        name: 'fno',
        label: 'FNO (neural operator, grid)',
        block: 'Grid',
        panelId: 't-fno-opts',
        fields: [
            { key: 'resolution', id: 't-fno-resolution', kind: 'triple', required: true },
            { key: 'squeeze', id: 't-fno-squeeze', kind: 'axis' },
            { key: 'latent_channels', id: 't-fno-latent', kind: 'number' },
            { key: 'num_fno_modes', id: 't-fno-modes', kind: 'number' },
            { key: 'num_fno_layers', id: 't-fno-layers', kind: 'number' },
        ],
    },
    {
        name: 'afno',
        label: 'AFNO (2-D operator, grid)',
        block: 'Grid',
        panelId: 't-afno-opts',
        fields: [
            { key: 'resolution', id: 't-afno-resolution', kind: 'triple', required: true },
            { key: 'squeeze', id: 't-afno-squeeze', kind: 'axis', required: true },
            { key: 'patch_size', id: 't-afno-patch', kind: 'pair', required: true },
            { key: 'embed_dim', id: 't-afno-embed', kind: 'number' },
            { key: 'depth', id: 't-afno-depth', kind: 'number' },
        ],
    },
    {
        name: 'deeponet',
        label: 'DeepONet (parameters in, field out)',
        block: 'Operator',
        panelId: 't-deeponet-opts',
        requires: 'deeponet',
        fields: [
            { key: 'parameters', id: 't-deep-params', kind: 'names', required: true },
            { key: 'trunk', id: 't-deep-trunk', kind: 'choice' },
            { key: 'trunk_count', id: 't-deep-trunk-count', kind: 'number' },
            { key: 'width', id: 't-deep-width', kind: 'number' },
            { key: 'branch_layers', id: 't-deep-branch-layers', kind: 'number' },
        ],
    },
];

export const FAMILY_NAMES: readonly string[] = FAMILIES.map((f) => f.name);

/** The family for a `Model.Name`; an unknown name is an error, never a
 * silent fallback to some other family's form. */
export function familyFor(name: string): Family {
    const family = FAMILIES.find((f) => f.name === name);
    if (!family) {
        throw new Error(`meshio++: unknown model family '${name}' (known: ${FAMILY_NAMES.join(', ')})`);
    }
    return family;
}

function positiveInts(text: string, count: number, what: string, example: string): number[] {
    const parts = text
        .split(',')
        .map((part) => part.trim())
        .filter((part) => part.length > 0)
        .map((part) => Number(part));
    if (parts.length !== count || parts.some((v) => !Number.isInteger(v) || v <= 0)) {
        throw new Error(`meshio++: ${what} must be ${count} positive integers, e.g. ${example}`);
    }
    return parts;
}

/**
 * A control's raw text -> the argument value, or `undefined` to omit the
 * key (a blank optional number/axis means "the server's default").
 */
export function fieldValue(field: FamilyField, raw: string): unknown {
    const text = raw.trim();
    if (text === '') {
        if (field.required) throw new Error(`meshio++: ${field.key} is required`);
        return undefined;
    }
    switch (field.kind) {
        case 'number': {
            const value = Number(text);
            if (!Number.isFinite(value)) throw new Error(`meshio++: ${field.key} must be a number`);
            return value;
        }
        case 'triple':
            return positiveInts(text, 3, field.key, '16,16,16');
        case 'pair':
            return positiveInts(text, 2, field.key, '8,8');
        case 'axis': {
            const axis = Number(text);
            if (![0, 1, 2].includes(axis)) throw new Error(`meshio++: ${field.key} must be a world axis 0, 1 or 2`);
            return axis;
        }
        case 'names': {
            const names = text
                .split(',')
                .map((part) => part.trim())
                .filter((part) => part.length > 0);
            if (!names.length) throw new Error(`meshio++: ${field.key} must name at least one Metadata key`);
            return names;
        }
        case 'choice':
            return text;
    }
}

/**
 * The family's arguments from its controls' raw values (`read(id)` returns
 * a control's text). Only the chosen family's keys are produced.
 */
export function familyArgs(family: Family, read: (id: string) => string): Record<string, unknown> {
    const args: Record<string, unknown> = {};
    for (const field of family.fields) {
        const value = fieldValue(field, read(field.id));
        if (value !== undefined) args[field.key] = value;
    }
    return args;
}
