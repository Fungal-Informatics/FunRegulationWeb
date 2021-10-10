import { MigrationDefinitionWithName } from "umzug/lib/migrationsList";
import { zap } from "../helpers/sql";

const migration: MigrationDefinitionWithName = {
	name: "01-boards",
	async up() {
		await zap.transaction.single`
            create table 
                boards
            ( 
                id uuid primary key,
                "ownerId" uuid not null references users(id),
                name text not null,
				type text,
                "createdAt" timestamptz not null,
                "updatedAt" timestamptz not null,
				"isArchived" boolean not null
            );
        `;
	},
	async down() {
		await zap.transaction.single`
            drop table boards;
        `;
	},
};

export default migration;
